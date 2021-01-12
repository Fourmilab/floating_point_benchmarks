
all:
	@echo "Please make for a specified language: ada algol60 algol68"
	@echo "    c c++ chapel cobol erlang fortran freebasic go haskell"
	@echo "    pascal pli prolog rust scala simula swift"

ada:	FORCE
	( cd ada ; rm -f fbench.o fbench ; gnatmake $(ADAOPTS) fbench.adb )

algol60:	FORCE
	( cd algol60; make clean; make )

algol68:	FORCE
	( cd algol68; make clean; make )

c:	FORCE
	( cd c; make clean; make )

c++:    FORCE
	( cd c++; make clean; make )

chapel: FORCE
	(cd chapel; make clean; make )

cobol:	FORCE
	( cd cobol; make clean; make )

erlang: FORCE
	( cd erlang; make clean; make )

fortran:	FORCE
	( cd fortran ; make clean; make )

freebasic:      FORCE
	( cd freebasic ; make clean ; make )

go:	FORCE
	( cd go; make clean; make )

haskell:	FORCE
	( cd haskell; make clean; make )

java:	FORCE
	( cd java; make clean; make )

modula2:        FORCE
	( cd modula2; make clean; make )

pascal: FORCE
	( cd pascal; make clean; make )

pli:    FORCE
	( cd pli; make clean; make )

prolog: FORCE
	( cd prolog; make clean; make )

rust:	FORCE
	( cd rust; make clean; make time )

scala:	FORCE
	( cd scala; make clean; make; make time )

simula:	FORCE
	( cd simula; make clean; make )

swift:	FORCE
	( cd swift; make clean; make )

dist:	FORCE
	make clean
	rm -f fbench.zip
	zip -r fbench.zip Makefile \
	    c \
	    ada \
	    algol60 \
	    algol68 \
	    atlast \
	    c++ \
	    chapel \
	    cobol \
	    erlang \
	    forth \
	    fortran \
	    freebasic \
	    go \
	    haskell \
	    java \
	    javascript \
	    julia \
	    lisp \
	    lua \
	    mathematica \
	    mbasic \
	    modula2 \
	    pascal \
	    perl \
	    php \
	    pli \
	    prolog \
	    python \
	    ruby \
	    rust \
	    scala \
	    scilab \
	    simula \
	    smalltalk \
	    swift \
	    visualbasic

clean:	FORCE
	rm -f core core.* *.bak *.zip
	( cd algol60; make clean )
	( cd algol68; make clean )
	( cd c; make clean )
	( cd c++; make clean )
	( cd chapel; make clean )
	( cd cobol; make clean )
	( cd erlang; make clean )
	( cd fortran; make clean )
	( cd freebasic; make clean )
	@bash -c "type go" >/dev/null 2>/dev/null ; if [ $$? -eq 0 ] ; then ( cd go; make clean ) ; fi
	( cd haskell; make clean )
	( cd java; make clean )
	( cd lua; make clean )
	( cd modula2 ; make clean )
	( cd pascal; make clean )
	( cd php; make clean )
	( cd pli; make clean )
	( cd prolog; make clean )
	@bash -c "type cargo" >/dev/null 2>/dev/null ; if [ $$? -eq 0 ] ; then ( cd rust; make clean ) ; fi
	( cd scala; make clean )
	( cd simula; make clean )
	( cd swift; make clean )


FORCE:
