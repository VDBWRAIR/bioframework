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
# 'I' = phred 40, '!' = phred 0, '5' = phred 20
${fastq_content}             @trimid\nNNNNNNNNNNAAAAAAAAAANNNNNNNNNN\n+\n5555555555IIIIIIIIII!!!!!!!!!!
${interleave_content}        ${fastq_content}\n${fastq_content}
${in_fastq} =                ${TESTOUTPUTDIR}/trim_me.fastq
${testout} =                 ${TESTOUTPUTDIR}/output
${trimmed_fq} =              @trimid\nAAAAAAAAAA\n+\nIIIIIIIIII\n
${out_fastq} =               ${trimmed_fq}${trimmed_fq}
${tool} =                    ${jip_path}/trim_reads.jip
${jipdb} =                   ${TESTOUTPUTDIR}/jip.db

*** Test Cases *** 
# Normal execution
trim_reads paired trims Ns
    Create File                     ${in_fastq}                 ${interleave_content}
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} ${in_fastq} -p -o ${testout}.fastq -t    expected_rc=1
    Should Contain                  ${stdout}             --trim-n is not supported at this time

trim_reads paired trims low quality
    Create File                     ${in_fastq}                 ${interleave_content}
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} ${in_fastq} -p -o ${testout}.fastq -q 25 30
    Verify File                     ${testout}.fastq            ${out_fastq}

trim_reads paired removes bases
    Create File                     ${in_fastq}                 ${interleave_content}
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} ${in_fastq} -p -o ${testout}.fastq -x 10 -10
    Verify File                     ${testout}.fastq            ${out_fastq}

trim_reads paired accepts unpaired input
    Create File                     ${in_fastq}                 ${fastq_content}
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} ${in_fastq} -p -o ${testout}.fastq -x 10 -10
    Verify File                     ${testout}.fastq            ${trimmed_fq}

trim_reads handles stdin/stdout pipe
    Ensure Tool Accepts Stdin And Outputs Stdout    ${tool} -x 10 -10    ${interleave_content}    ${out_fastq}
# End Normal execution

# Ensure normalized usage statements
trim_reads supports input and output
    Normalized input And Output Usage Statement    ${tool}

# Expected failures
# End Expected failures
