#!/usr/bin/env cwl-runner

cwlVersion: "cwl:draft-3.dev3"

requirements:
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      - envName: PATH
        envValue: 'myconda/bin:\$PATH'

class: CommandLineTool

inputs:
  - id: apps
    type:
      type: array
      items: string
    inputBinding: { position: 1 }

outputs:
  - id: output
    type: File
    outputBinding:
        outputEval: |
            ${
                var pkgPaths = [];
                for(var input in inputs) {
                    pkgPaths.push('myconda/bin' + input)
                }
                return pkgPaths;
            }

baseCommand: [conda, install]
