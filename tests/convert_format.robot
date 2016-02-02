*** Settings ***
Library             Process
Library             OperatingSystem
Library             Collections 
Library             String
Suite Teardown      Terminate All Processes    

*** Variables ***
${jip_dir} =                 ${CURDIR}/../pyjip
${in_fastq} =                ${CURDIR}/testinput/test.fastq
${testout} =                 ${CURDIR}/testoutput/output
${fastaout} =                >id1\nATGC
${tool} =                    ${jip_dir}/convert_format.jip

*** Test Cases *** 
convert_format file to stdout
    Append To Environment Variable  PYTHONPATH    ${CURDIR}/..
    ${result} =    Run Process      ${tool} -i ${in_fastq} --out-format fasta     shell=True    stderr=STDOUT
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}            0 
    Should Be Equal As Strings      ${result.stdout}        ${fastaout}

convert_format file to file
    ${result} =    Run Process      rm ${testout}.fasta; ${tool} -i ${in_fastq} -o ${testout}.fasta    shell=True    stderr=STDOUT
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}            0 
    ${actual_contents} =            Get File                ${testout}.fasta
    ${actual_contents} =            Strip String            ${actual_contents}
    Should Be Equal As Strings      ${fastaout}             ${actual_contents}

convert_format stdin to stdout
    ${result} =    Run Process      cat ${in_fastq} | ${tool} --in-format fastq --out-format fasta     shell=True    stderr=STDOUT
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}            0 
    Should Be Equal As Strings      ${result.stdout}    ${fastaout}

convert_format missing outformat
    ${result} =    Run Process      ${tool} -i ${in_fastq}    shell=True    stderr=STDOUT
    Log To Console                  ${result.stdout}
    ${result} =    Run Process      ${tool} -i ${in_fastq}    shell=True    stderr=STDOUT
    Should Be Equal As Integers     ${result.rc}            1 
    ${lines}                        Get Regexp Matches      ${result.stdout}        Have to supply output format if output is stdout
    Should Be Equal As Strings      @{lines}[0]             Have to supply output format if output is stdout

convert_format missing informat
    ${result} =    Run Process      cat ${in_fastq} | ${tool} --out-format fasta    shell=True    stderr=STDOUT
    Log To Console                  ${result.stdout}
    Should Be Equal As Integers     ${result.rc}            1 
    ${lines}                        Get Regexp Matches      ${result.stdout}        Have to supply input format if input is stdin
    Should Be Equal As Strings      @{lines}[0]             Have to supply input format if input is stdin
