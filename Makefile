#!/usr/bin/make -f

probeAnnotation21kdatMethUsed.csv:
	wget -O $@ "http://labs.genetics.ucla.edu/horvath/dnamage/$@"
