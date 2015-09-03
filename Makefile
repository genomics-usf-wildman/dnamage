#!/usr/bin/make -f

all: probeAnnotation21kdatMethUsed.csv AdditionalFile3.csv datMiniAnnotation27k.csv

probeAnnotation21kdatMethUsed.csv:
	wget -O $@ "http://labs.genetics.ucla.edu/horvath/dnamage/$@"


AdditionalFile3.csv:
	wget -O $@ "http://labs.genetics.ucla.edu/horvath/dnamage/$@"

datMiniAnnotation27k.csv:
	wget -O $@ "http://labs.genetics.ucla.edu/horvath/dnamage/$@"
