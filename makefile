
JFLAGS = -g
JC = javac
JAR = jar
JVM = 1.8
out_jarfile = tfbs_scan.jar

BIN = classes

LIB = lib/commons-cli-1.3.1.jar

JAVAFLAGS = -g -d $(BIN) -cp $(LIB) -target $(JVM)

SRC = src/osu/megraw/llscan/*.java

manifest_file = MANIFEST.MF

compile:
	$(JC) $(JAVAFLAGS) $(SRC)

makejar: 
	$(JAR) cvfm $(out_jarfile) $(manifest_file) -C $(BIN) .

all: compile makejar
	echo "$(out_jarfile) created!"

clean:
	rm -f -r $(BIN)/*
	rm -f $(out_jarfile)
