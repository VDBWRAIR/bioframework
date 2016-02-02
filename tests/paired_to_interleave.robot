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
    ${result} =    Run Process      ${tool} -f ${in_fastq} -r ${in_fastq}     shell=True    stderr=STDOUT
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}            0 
    ${out_no_newline}               Strip String            ${out_fastq}${out_fastq}
    Should Be Equal As Strings      ${result.stdout}        ${out_no_newline}

paired_to_interleave file to file
    ${result} =    Run Process      ${tool} -f ${in_fastq} -r ${in_fastq} -o ${testout}.fastq     shell=True    stderr=STDOUT
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}            0 
    ${actual_contents} =            Get File                ${testout}.fastq
    Should Be Equal As Strings      ${out_fastq}${out_fastq}            ${actual_contents}

paired_to_interleave outformat specified
    ${result} =    Run Process      ${tool} -f ${in_fastq} -r ${in_fastq} -o ${testout}.fasta --out-format fasta    shell=True    stderr=STDOUT
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}                0 
    ${actual_contents} =            Get File                    ${testout}.fasta
    Should Be Equal As Strings      ${out_fasta}                ${actual_contents}
# End Normal execution

# Expected failures
paired_to_interleave outformat not fastq or fasta raises error
    ${result} =    Run Process      ${tool} -f ${in_fastq} -r ${in_fastq} --out-format foo    shell=True    stderr=STDOUT
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}            1 
    ${lines}                        Get Regexp Matches      ${result.stdout}        output format can only be fasta or fastq. You provided 'foo'
    Should Be Equal As Strings      @{lines}[0]             output format can only be fasta or fastq. You provided 'foo'
# End Expected failures
