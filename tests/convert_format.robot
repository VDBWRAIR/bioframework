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
${fastaout} =                >id1\nATGC
${tool} =                    convert_format.jip

*** Test Cases *** 
convert_format file to stdout
    Append To Environment Variable  PYTHONPATH    ${CURDIR}/..
    ${result} =    Run Shell      ${tool} -i ${in_fastq} --out-format fasta 
    Should Be Equal As Strings      ${result.stdout}        ${fastaout}

convert_format file to file
    ${result} =    Run Shell      ${tool} -i ${in_fastq} --out-format fasta   
    Should Be Equal As Strings      ${result.stdout}        ${fastaout}

convert_format file to file
    ${result} =    Run Shell      ${tool} -i ${in_fastq} -o ${testout}.fasta 
    ${actual_contents} =            Get File                ${testout}.fasta
    ${actual_contents} =            Strip String            ${actual_contents}
    Should Be Equal As Strings      ${fastaout}             ${actual_contents}

convert_format stdin to stdout
    ${result} =    Run Shell      cat ${in_fastq} | ${tool} --in-format fastq --out-format fasta 
    Should Be Equal As Strings      ${result.stdout}    ${fastaout}

convert_format missing outformat
    ${result} =    Command Should Fail      ${tool} -i ${in_fastq}    
    Should Be SubString         ${result.stdout}        Have to supply output format if output is stdout

convert_format missing informat
    ${result} =    Command Should Fail      cat ${in_fastq} | ${tool} --out-format fasta
        Should Be Substring     ${result.stdout}        Have to supply input format if input is stdin 
