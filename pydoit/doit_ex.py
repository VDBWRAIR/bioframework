from fn import _
from helpers import D, F, fifo, to_fifo, group_fastqs, swap_ext
from glob import glob1
#NOTE:
'''
 they're dropping python2 support in next release
 it's by modification time
 the output will always have to be created as a fifo, but not the input
 to wrap "targets" in fifo
 'task_dep' does exist
'''

'''Options'''
class Opts(object): pass
opts = Opts()
opts.cutadapt = { 'headcrop' : 33 }
opts.bwa = { 'threads' : 1, 'keep_temp' : False}
opts.tagreads = { 'CA' : "Foo"}
opts.freebayes = { 'observations' : 5, 'haplotype_length' : 50}

'''Input files'''
DIR="."
opts.ref = F("ref.fasta")
R1s, R2s, unpaireds =  group_fastqs(DIR)
R1, R2, unpaired, paired, merged, consensus = map(F, ["R1", "R2", "unpaired", "paired", "merged", "consensus"])
sffs = glob1(DIR, "*.sff")
unpaireds += [swap_ext(sff, 'fastq') for sff in sffs]

'''example tasks'''
def task_index_ref():
    return { 'targets' : [opts.ref.fai],
             'file_dep' : [opts.ref],
             'actions' : ['bwa index %(dependencies)'] }

def task_sff2fastq(): 
    @fifo('o')
    def convert_sff(i, o): SeqIO.convert(i, 'sff', o, 'fastq')

    for sff in sffs:
        yield {'target' : [swap_ext(sff, 'fastq')],
               'file_deps' : [sff],
               'actions' : convert_sff}

''' examples of programmatically created tasks 
tasks must be wrapped in functions starting with `task_`
(so using lambda: dict) '''


def cutadapt_paired(dependencies, targets, quality):
    sh.cutadapt(o=targets[0], p=targets[1], q="{0},{0}".format(quality), *dependencies)

def cutadapt_up(dependencies, targets, quality):
    sh.cutadapt(q=quality, o=targets[0], *dependencies) 
    #'actions': ["cutadapt %(dependencies)"] # not work

def gen_cutdapt(fqs):
   func = cutadapt_paired if len(fqs) > 1 else cutadapt_up
   return { 'targets'  : map(_.cutadapt, fqs),
            'file_dep' : fqs, 
            'actions' : [to_fifo(func, 'targets'), [], opts.cutadapt] }

def task_run_cutadapt():
    for fqs in zip(R1s, R2s) + unpaireds:
        yield gen_cutdapt(fqs) 


merge_files = D({
    'actions' : ['cat %(dependencies) > %(targets)']})

# _.cutdadapt = lambda x: x.cudadapt
# `.assoc` is like `.update` but returns a new dictionary without changing the old one
task_R1 = lambda: merge_files.assoc(file_dep=map(_.cutadapt, R1s), targets=R1.fastq)
task_R2 = lambda: merge_files.assoc(file_dep=map(_.cutadapt, R2s), targets=R2.fastq)
task_unpaired  = lambda: merge_files.assoc(file_dep=map(_.cutadapt, unpaireds), targets=unpaired.fastq)


@fifo('targets')
def bwa_map(dependencies, targets, threads, keep_temp):
   '''`dependencies` and `targets` are lists auto-provided by
   pydoit with `file_deps` and `targets` from the task dictionary.
   The keyword arguments are provided by the dictionary in the
   third element of `'actions'`, in this case `opts.bwa`'''
   sh.bwa.mem(dependencies, t=threads, _out=targets[0])

mapping =   D({'file_dep' : [opts.ref],
             'actions' : [bwa_map, [], opts.bwa]}) #opts.bwa here sends the kwargs to `bwa_map`

#`.apply(key=func)` applies the func to the value at that key

mapping_up = mapping.apply(file_dep = _ + [unpaired.fastq]).assoc(targets=[unpaired.bam])
mapping_paired = mapping.apply(file_dep = _ + [R1.fastq, R2.fastq]).assoc(targets=[paried.bam])
task_maping_up = lambda : mapping_up
task_maping_paired = lambda : mapping_paired 

def task_merge_bam():
    return { 'targets' : [unsorted.bam]
             'file_dep' : [paired.bam, unpaired.bam], #'task_dep' : 'mapping',
             'actions' : ['mkfifo %(targets)', 'samtools merge %(targets) %(dependencies)']}

def make_task(target, dep, action):
    return { 'targets' : [target], 'file_dep' : [dep], 'actions' : [action] }

task_sort_bam = lambda : make_task(merged.bam, unsorted.bam, 'mkfifo %(targets) && samtools sort %(dependencies) > %(targets)')
task_index_bam = lambda : make_task(merged.bam.bai, merged.bam, 'samtools index %(dependencies)')
task_tag_bam = lambda : { 'file_dep' : [merged.bam],
                          'actions' : ['tagreads %(dependencies) --CN %s' % opts.tagreads['CN'] ] }
# this is screwed up because the tagged bam file may not be sorted. Should swap out `tagreads` for 
# something that creates a new file, it's simpler
# also that doesn't play well with fifo

@fifo('targets')
def run_freebayes(dependencies, targets, haplotype_length, observations):
    sh.freebayes(dependencies[0], f=dependencies[1], 
            haplotype_length=haplotype_length, C=observations, _out=targets[0])

def task_freebayes(): 
    return {  'targets' : [merged.bam.vcf]
              'file_dep' : [merged.bam, opts.ref, ref.fasta.fai]
              'task_dep' : ['tag_bam']
              'actions' : [run_freebayes, [], opts.freebayes] }

def task_consensus():
    def TODO(dependecies, targets):
        pass
    return { 'targets' : [consensus.fasta]
             'file_dep' : [merged.bam.vcf, opts.ref]
             'actions' : [TODO] }

'''
actions take kwargs. so add the dict associated with a task (e.g., ngs_filter)
dynamically, based on func name #
dy of task-creators are executed even if the task is not going to be executed.
'''
