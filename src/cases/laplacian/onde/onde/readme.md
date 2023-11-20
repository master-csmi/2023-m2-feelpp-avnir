# Run the test case Onde

Click on 'build' with preset 'default'. Then run, when located at the root of the project, the following command:

```bash
 cd build/default/src
./feelpp_2mfs_laplacian --config-file ../../../src/cases/laplacian/onde/onde/onde.cfg
```

## Notes

### json changes

Some tests needs changes in parameters, these lines of `onde.json` :
```json
        "Omega":{
            "c": "4",
            "x0":"1",
            "y0":"1",
            "sigma":"0.05",
            "a": "0.3"
        }
```

Became:
```json
        "Omega":{
            "c": "4"
        }
```