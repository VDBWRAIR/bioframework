from fn import _
from helpers import D, F, fifo, to_fifo, group_fastqs
#NOTE:
'''
 they're dropping python2 support in next release
 it's by modification time
 the output will always have to be created as a fifo, but not the input
 to wrap "targets" in fifo
 'task_dep' does exist
'''

class Opts(object): pass
opts = Opts()
opts.cutadapt = { 'headcrop' : 33 }
opts.bwa = { 'threads' : 1, 'keep_temp' : False}
opts.ref = F("ref.fasta")

'''example tasks'''
def task_index_ref():
    return { 'file_dep' : [opts.ref],
             'actions' : ['bwa index %(dependencies)'],
             'targets' : [opts.ref.fai] }
def task_merge_bam():
    return { 'task_dep' : 'mapping'}


''' examples of programmatically created tasks '''
merge_files = D({
    'actions' : ['cat %(dependencies) %(targets)']})

# _.cutdadapt = lambda x: x.cudadapt
task_R1 = lambda: merge_files.assoc(file_dep=map(_.cutadapt, R1s), targets=R1.fastq)
task_R2 = lambda: merge_files.assoc(file_dep=map(_.cutadapt, R2s), targets=R2.fastq)
task_unpaired  = lambda: merge_files.assoc(file_dep=map(_.cutadapt, unpaireds), targets=unpaired.fastq)

DIR="."
R1s, R2s, unpaireds =  group_fastqs(DIR)
R1, R2, unpaired, paired = map(F, ["R1", "R2", "unpaired", "paired"])


mapping =   D({'file_dep' : [opts.ref],
             'actions' : [bwa_map, [], opts.bwa]})

mapping_up = mapping.apply(file_dep = _ + [unpaired.fastq]).assoc(targets=[unpaired.bam])
mapping_paired = mapping.apply(file_dep = _ + [R1.fastq, R2.fastq]).assoc(targets=[paried.bam])
task_maping_up = lambda : mapping_up
task_maping_paired = lambda : mapping_paired



@fifo('targets')
def bwa_map(dependencies, targets, threads, keep_temp):
   sh.bwa.mem(*dependencies, t=threads, _out=targets[0])

R1s, R2s, unpaireds = group_fastqs(files)

def task_run_cutadapt():
    for fqs in zip(R1s, R2s) + unpaireds:
        yield gen_cutdapt(fqs)

def gen_cutdapt(*fqs):
   func = cutadapt_paired if len(fqs) > 1 else cutadapt_up
   return {
    'file_dep' : fqs,
   'targets'  : map(_.cutadapt, fqs),
   'actions' : [to_fifo(func, 'targets'), [], opts.cutadapt]}

def cutadapt_paired(dependencies, targets, quality):
    sh.cutadapt(o=targets[0], p=targets[1], q="{0},{0}".format(quality), *dependencies)

def cutadapt_up(dependencies, targets, quality):
    sh.cutadapt(q=quality, o=targets[0], *dependencies)

    #'actions': ["cutadapt %(dependencies)"] # not work

#fwd, rev <- fwd, rev : cutadapt_paired


def task_sff2fastq():
    convert_sff = lambda i,o: SeqIO.convert(i, 'sff', o, 'fastq')
    convert_sff = to_fifo(convert_sff, 'o')
    for sff in sffs:
        yield {'target' : [swap_ext(sff, 'fastq')],
               'file_deps' : [sff],
               'actions' : convert_sff}
'''
actions take kwargs. so add the dict associated with a task (e.g., ngs_filter)
dynamically, based on func name #
dy of task-creators are executed even if the task is not going to be executed.
'''
