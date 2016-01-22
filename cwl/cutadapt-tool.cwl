#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

requirements:
    - class: EnvVarRequirement
      envDef:
        - envName: PATH
          envValue: "myconda/bin:$PATH"

inputs:
    - id: reads
      type:
        type: array
        items: File
      inputBinding: { position: 3 }
      description: "Input fastq files"

    - id: qualcutoff
      type:
        type: array
        items: string
      inputBinding:
        itemSeparator: ","
        prefix: -q
      default:
        - '25'
        - '25'

outputs:
    - id: out1
      type: File
      outputBinding: { glob: out1.cutadapt }

    - id: out2
      type: File
      outputBinding: { glob: out2.cutadapt }

baseCommand: [myconda/bin/cutadapt]

arguments:
    - prefix: -o 
      position: 1
      valueFrom: out1.cutadapt

    - prefix: -p 
      position: 2
      valueFrom: out2.cutadapt

stdout: output.sam
