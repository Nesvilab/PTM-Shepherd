cd ../../
./gradlew packageWDeps
cd test/shepherd/
java -jar ../../build/libs/ptmshepherd-0.2.7.jar config_test.txt 

