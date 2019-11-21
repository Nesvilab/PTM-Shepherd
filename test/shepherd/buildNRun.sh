cd ../../
./gradlew packageWDeps
cd test/shepherd/
java -jar ../../build/libs/ptmshepherd-0.2.14.jar shepherd_smooth5.config 


