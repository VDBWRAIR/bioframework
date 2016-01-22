#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: "cwl:draft-3.dev3"

inputs:
    - id: reads
      type:
        type: array
        items: File

outputs:
    - id: read1.cutadapt
      type: "null"

requirements:
    - class: SubworkflowFeatureRequirement
    - class: EnvVarRequirement
      envDef:
        - envName: PATH
          envValue: "myconda/bin:$PATH"

steps:
    - id: trimming
      run: cutadapt-tool.cwl
      inputs:
        - { id: reads, source: "#reads" }
      outputs:
        - { id: out1 }
        - { id: out2 }
