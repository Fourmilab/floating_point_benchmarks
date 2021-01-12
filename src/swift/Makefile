
BUILDOPTS = -O

fbench: fbench.swift
	swiftc $(BUILDOPTS) fbench.swift
clean:
	rm -f fbench

#	Time compiled code

time:	fbench
	time -p ./fbench

#	Time in interpretive execution environment

timeREPL:
	time -p swift fbench.swift

