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
${trimmed_fq} =              @trimid\nAAAAAAAAAA\n+\nIIIIIIIIII
${out_fastq} =               ${trimmed_fq}\n${trimmed_fq}\n
${tool} =                    ${jip_path}/trim_reads.jip
${jipdb} =                   ${TESTOUTPUTDIR}/jip.db

*** Test Cases *** 
# Normal execution
trim_reads paired trims Ns
    Create File                     ${in_fastq}                 ${interleave_content}
    ${result} =    Run Process      ${tool} ${in_fastq} -p -o ${testout}.fastq -t     shell=True    stderr=STDOUT    env:JIP_DB=${jipdb}
    Log To Console                  ${result.stdout}
    Should Contain                  ${result.stdout}             --trim-n is not supported at this time
    #Should Be Equal As Integers     ${result.rc}                0 
    #Verify File                     ${testout}.fastq            ${out_fastq}

trim_reads paired trims low quality
    Create File                     ${in_fastq}                 ${interleave_content}
    ${result} =    Run Process      ${tool} ${in_fastq} -p -o ${testout}.fastq -q 25 30    shell=True    stderr=STDOUT    env:JIP_DB=${jipdb}
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}                0 
    Verify File                     ${testout}.fastq            ${out_fastq}

trim_reads paired removes bases
    Create File                     ${in_fastq}                 ${interleave_content}
    ${result} =    Run Process      ${tool} ${in_fastq} -p -o ${testout}.fastq -x 10 -10    shell=True    stderr=STDOUT    env:JIP_DB=${jipdb}
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}                0 
    Verify File                     ${testout}.fastq            ${out_fastq}

trim_reads paired accepts interleave input
    Create File                     ${in_fastq}                 ${interleave_content}
    ${result} =    Run Process      ${tool} ${in_fastq} -p -o ${testout}.fastq -x 10 -10    shell=True    stderr=STDOUT    env:JIP_DB=${jipdb}
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}                0 
    Verify File                     ${testout}.fastq            ${out_fastq}

trim_reads handles stdin/stdout pipe
    Create File                     ${in_fastq}                 ${interleave_content}
    ${result} =    Run Process      cat ${in_fastq} | ${tool} -x 10 -10 | cat > ${testout}.fastq  shell=True    stderr=STDOUT    env:JIP_DB=${jipdb}
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}                0 
    Verify File                     ${testout}.fastq            ${out_fastq}

# End Normal execution

# Expected failures
# End Expected failures
