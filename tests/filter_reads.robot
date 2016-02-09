*** Settings ***
Library             Process
Library             OperatingSystem
Library             Collections 
Library             String
Suite Teardown      Common Teardown
Test Setup          Common Setup
Resource            common_robot.txt

*** Variables ***
${jip_path} =                jip_modules
# 'I' = phred 40, '!' = phred 0, '5' = phred 20
${testout1} =                 ${TESTOUTPUTDIR}/output
${testout2} =                 ${TESTOUTPUTDIR}/output2
${in_fastq} =       tests/testinput/filter_me.fastq
${expected_fastq} =     tests/testinput/expected__R2__.fastq
${index1} =     tests/testinput/filter_me_I1_.fastq
${index2} =     tests/testinput/filter_me_I2_.fastq 
${dropNs} =     ${jip_path}/drop_ns.jip
${filterIndex} =        ${jip_path}/interleaved_index_filter.jip
${cmd} =        ${filterIndex} --minimum 32 --index1 ${index1} --index2 ${index2} | ${dropNs} ${in_fastq} -p --max-n 0 

*** Test Cases *** 
# Normal execution
40 Quality doesn't change file
        Run JIP Tool And Return RC And Output       ${filterIndex} ${in_fastq} -o ${testout2}.fastq --minimum 0 --index1 ${index1} --index2 ${index2} 
        Files Should Be Equal   ${testout2}.fastq       ${in_fastq}

High maximum N's doesn't change file
    Run JIP Tool And Return RC And Output       ${dropNs} ${in_fastq} -p -o ${testout1}.fastq  --max-n 999
    Files Should Be Equal   ${testout1}.fastq       ${in_fastq}

Works on expected file
    Run JIP Tool And Return RC And Output      ${dropNs} ${in_fastq} -p -o ${testout1}.fastq --max-n 0
    Run JIP Tool And Return RC And Output       ${filterIndex} ${testout1}.fastq -o ${testout2}.fastq --minimum 32 --index1 ${index1} --index2 ${index2} -o ${testout2}
        Files Should Be Equal   ${testout2}     ${expected_fastq} 

accepts stdin and stdout
        ${input} =      Get File        ${in_fastq}
        ${expected} =      Get File        ${expected_fastq}
        Ensure Tool Accepts Stdin And Outputs Stdout    ${cmd}  ${input}        ${expected}


# End Normal execution 
# Expected failures
# End Expected failures
