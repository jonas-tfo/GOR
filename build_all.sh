rm -rf out/
javac --release 21 -d out $(find src/main -name "*.java")
jar cfm train.jar train_manifest.txt -C out . -C src/main .
jar cfm predict.jar predict_manifest.txt -C out . -C src/main .
jar cfm validate.jar validate_manifest.txt -C out . -C src/main .
tar -czvf tar.gz train.jar predict.jar validate.jar
mkdir -p out/artifacts
mv train.jar out/artifacts/
mv predict.jar out/artifacts/
mv validate.jar out/artifacts/
mv tar.gz out/artifacts/
