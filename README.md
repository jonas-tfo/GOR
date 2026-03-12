# GOR

An implementation of the GOR (Garnier-Osguthorpe-Robson) protein secondary structure prediction algorithm in Java, supporting methods GOR1 through GOR5. Predicts three secondary structure states — coil (C), beta-strand (E), and helix (H) for given amino acid sequences.

## Build

Requires Java 21.

```bash
./build_all.sh
```

Produces three JARs in `out/artifacts/`, for training, prediction and validation respectively.

## Usage

**Training**
```bash
java -jar train.jar --db <training_file> --method <gor1|gor3|gor4|gor5> --model <output_model>
```

**Prediction**
```bash
java -jar predict.jar --model <model_file> --seq <sequence_file> [--probabilities] [--colored]
```

**Validation**
```bash
java -jar validate.jar --db <validation_file> --method <gor1|gor3|gor4|gor5> --model <model_file>
```

---

## License

MIT License

Copyright (c) 2026

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
