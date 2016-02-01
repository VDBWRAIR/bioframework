import jsonschema
import json
example = json.loads(open('example.json'))
schema = json.loads(open('schema.json'))
jsonschema.validate(example, schema)
