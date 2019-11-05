cd ../../
./gradlew packageWDeps
cd test/shepherd/
java -jar ../../build/libs/ptmshepherd-0.2.12.jar shepherd_smooth5.config 

echo ../../../PTMShepherdTestsFromAlexey/FromAlexey/Aug9_open_v2/shepherd.config 

echo shepherd_smooth5.config

