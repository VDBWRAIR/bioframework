#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3.dev3"

class: CommandLineTool

inputs:
  - id: reference
    type: File
    inputBinding: { position: 2 }
    description: "Reference file path"

  - id: reads
    type:
      type: array
      items: File
    inputBinding: { position: 3 }
    description: "Input fastq files"

outputs:
  - id: sam
    type: File
    outputBinding: { glob: output.sam }

baseCommand: [bwa, mem]

arguments:
  - valueFrom: $(runtime.cores)
    position: 1
    prefix: -t

stdout: output.sam
