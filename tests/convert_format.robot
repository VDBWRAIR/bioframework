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
${tool} =                    ${jip_path}/convert_format.jip

*** Test Cases *** 
convert_format file to stdout
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} -i ${in_fastq} --out-format fasta
    Should Be Equal As Strings      ${stdout}        ${fastaout}

convert_format file to file
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} -i ${in_fastq} -o ${testout}.fasta
    Verify File    ${testout}.fasta    ${fastaout}

convert_format stdin to stdout
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    cat ${in_fastq} | ${tool} --in-format fastq --out-format fasta
    Should Be Equal As Strings      ${stdout}    ${fastaout}

convert_format missing outformat
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    ${tool} -i ${in_fastq}     expected_rc=1
    ${lines}                        Get Regexp Matches      ${stdout}        Have to supply output format if output is stdout
    Should Be Equal As Strings      @{lines}[0]             Have to supply output format if output is stdout

convert_format missing informat
    ${rc}    ${stdout} =            Run JIP Tool And Return RC And Output    cat ${in_fastq} | ${tool} --out-format fasta    expected_rc=1
    ${lines}                        Get Regexp Matches      ${stdout}        Have to supply input format if input is stdin
    Should Be Equal As Strings      @{lines}[0]             Have to supply input format if input is stdin
