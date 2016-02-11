*** Settings ***
Library             Process
Library             OperatingSystem
Library             Collections 
Library             String
Suite Teardown      Common Teardown
Test Setup          Common Setup
Resource            common_robot.txt

*** Variables ***
${jip_path} =                ${PROJDIR}/jip_modules
${in_fastq} =                ${TESTINPUTDIR}/test.fastq
${testout} =                 ${TESTOUTPUTDIR}/output
${out_fastq} =               @id1\nATGC\n+\n!!!!\n
${out_fasta} =               >id1\nATGC\n>id1\nATGC\n
${tool} =                    ${jip_path}/paired_to_interleave.jip

*** Test Cases *** 
# Normal execution
paired_to_interleave file to stdout
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} ${in_fastq} ${in_fastq}
    Compare Stripped Contents       ${out_fastq}${out_fastq}    ${stdout}

paired_to_interleave file to file
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} ${in_fastq} ${in_fastq} -o ${testout}.fastq
    Verify File     ${testout}.fastq    ${out_fastq}${out_fastq}

paired_to_interleave outformat specified
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} ${in_fastq} ${in_fastq} -o ${testout}.fasta --out-format fasta
    Verify File    ${testout}.fasta    ${out_fasta}
# End Normal execution

# Ensure normalized usage statements
trim_reads supports input and output
    Normalized input And Output Usage Statement    ${tool}    supports_paired=False
# End Normalized usage

# Expected failures
paired_to_interleave outformat not fastq or fasta raises error
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} ${in_fastq} ${in_fastq} --out-format foo    expected_rc=1
    ${lines}                        Get Regexp Matches      ${stdout}        output format can only be fasta or fastq. You provided 'foo'
    Should Be Equal As Strings      @{lines}[0]             output format can only be fasta or fastq. You provided 'foo'

paired_to_interleave must be given two input files
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} ${in_fastq} --out-format fastq    expected_rc=1
    Should Contain    ${stdout}     Must supply two input files to convert to interleave

# End Expected failures
