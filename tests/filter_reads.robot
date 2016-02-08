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
${in_fastq} =                testinput/filter_me.fastq
${testout1} =                 ${TESTOUTPUTDIR}/output
${testout2} =                 ${TESTOUTPUTDIR}/output2
${expected_fastq} =     testinput/expected__R2__.fastq
${jipdb} =                   ${TESTOUTPUTDIR}/jip.db
${index1} =     testinput/filter_me_I1_.fastq
${index2} =     testinput/filter_me_I2_.fastq 
${dropNs} =     drop_ns.jip
${filterIndex} =        interleaved_index_filter.jip
*** Test Cases *** 
# Normal execution
40 Quality doesn't change file
        Run JIP Tool And Return RC And Output       ${filterIndex} ${in_fastq} -o ${testout2}.fastq --minimum 0 --index1 ${index1} --index2 ${index2} -o ${testout2}
        Files Should Be Equal   ${testout2}.fastq       ${in_fastq}

High maximum N's doesn't change file
    Run JIP Tool And Return RC And Output
      ${dropNs} ${in_fastq} -p -o ${testout1}.fastq -t --max-n 999
    Files Should Be Equal   ${testout1}.fastq       ${in_fastq}

filter out Ns
    Run JIP Tool And Return RC And Output      ${dropNs} ${in_fastq} -p -o ${testout1}.fastq -t --max-n 1
    Run JIP Tool And Return RC And Output       ${filterIndex} ${testout1}.fastq -o ${testout2}.fastq --minimum 32 --index1 ${index1} --index2 ${index2} -o ${testout2}
        Files Should Be Equal   ${testout2}     ${expected_fastq} 
# End Normal execution

# Expected failures
# End Expected failures
