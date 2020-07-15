# Makefile for SWMM5

objs = swmm5.o climate.o controls.o culvert.o datetime.o dwflow.o dynwave.o error.o \
       exfil.o findroot.o flowrout.o forcmain.o gage.o gwater.o hash.o hotstart.o iface.o \
       infil.o inflow.o input.o inputrpt.o keywords.o kinwave.o landuse.o lid.o \
       lidproc.o link.o massbal.o mathexpr.o mempool.o node.o odesolve.o output.o \
       project.o qualrout.o rain.o rdii.o report.o roadway.o routing.o runoff.o shape.o snow.o \
       stats.o statsrpt.o subcatch.o surfqual.o table.o toposort.o transect.o treatmnt.o xsect.o


swmm5 : $(objs)
	cc -o libswmm5.so $(objs) -fopenmp -lm -lpthread -shared

swmm5.o       : consts.h macros.h enums.h error.h datetime.h objects.h funcs.h text.h globals.h swmm5.h
	gcc -c -fPIC swmm5.c
climate.o     : headers.h
	gcc -c -fPIC climate.c
controls.o    : headers.h
	gcc -c -fPIC controls.c
culvert.o     : headers.h
	gcc -c -fPIC culvert.c
datetime.o    : datetime.h
	gcc -c -fPIC datetime.c
dwflow.o      : headers.h
	gcc -c -fPIC dwflow.c
dynwave.o     : headers.h
	gcc -c -fPIC dynwave.c
error.o       : error.h
	gcc -c -fPIC error.c
exfil.o       : headers.h infil.h exfil.h
	gcc -c -fPIC exfil.c
findroot.o    : findroot.h
	gcc -c -fPIC findroot.c
flowrout.o    : headers.h
	gcc -c -fPIC flowrout.c
forcmain.o    : headers.h
	gcc -c -fPIC forcmain.c
gage.o        : headers.h
	gcc -c -fPIC gage.c
gwater.o      : headers.h odesolve.h
	gcc -c -fPIC gwater.c
hash.o        : hash.h
	gcc -c -fPIC hash.c
hotstart.o      : headers.h
	gcc -c -fPIC hotstart.c
iface.o       : headers.h
	gcc -c -fPIC iface.c
infil.o       : headers.h infil.h
	gcc -c -fPIC infil.c
inflow.o      : headers.h
	gcc -c -fPIC inflow.c
input.o       : headers.h lid.h
	gcc -c -fPIC input.c
inputrpt.o    : headers.h lid.h
	gcc -c -fPIC inputrpt.c
keywords.o    : text.h
	gcc -c -fPIC keywords.c
kinwave.o     : headers.h findroot.h
	gcc -c -fPIC kinwave.c
landuse.o     : headers.h
	gcc -c -fPIC landuse.c
lid.o         : headers.h infil.h lid.h
	gcc -c -fPIC lid.c
lidproc.o     : headers.h lid.h
	gcc -c -fPIC lidproc.c
link.o        : headers.h
	gcc -c -fPIC link.c
massbal.o     : headers.h
	gcc -c -fPIC massbal.c
mathexpr.o    : mathexpr.h
	gcc -c -fPIC mathexpr.c
mempool.o     : mempool.h
	gcc -c -fPIC mempool.c
node.o        : headers.h findroot.h
	gcc -c -fPIC node.c
odesolve.o    : odesolve.h
	gcc -c -fPIC odesolve.c
output.o      : headers.h
	gcc -c -fPIC output.c
project.o     : headers.h hash.h lid.h mempool.h
	gcc -c -fPIC project.c
qualrout.o    : headers.h
	gcc -c -fPIC qualrout.c
rain.o        : headers.h
	gcc -c -fPIC rain.c
rdii.o        : headers.h
	gcc -c -fPIC rdii.c
report.o      : headers.h
	gcc -c -fPIC report.c
roadway.o     : headers.h
	gcc -c -fPIC roadway.c
routing.o     : headers.h
	gcc -c -fPIC routing.c
runoff.o      : headers.h odesolve.h
	gcc -c -fPIC runoff.c
shape.o       : headers.h
	gcc -c -fPIC shape.c
snow.o        : headers.h
	gcc -c -fPIC snow.c
stats.o       : headers.h
	gcc -c -fPIC stats.c
statsrpt.o    : headers.h lid.h
	gcc -c -fPIC statsrpt.c
subcatch.o    : headers.h lid.h odesolve.h
	gcc -c -fPIC subcatch.c
surfqual.o    : headers.h lid.h
	gcc -c -fPIC surfqual.c
table.o       : headers.h
	gcc -c -fPIC table.c
toposort.o    : headers.h
	gcc -c -fPIC toposort.c
transect.o    : headers.h
	gcc -c -fPIC transect.c
treatmnt.o    : headers.h
	gcc -c -fPIC treatmnt.c
xsect.o       : headers.h findroot.h
	gcc -c -fPIC xsect.c

clean:
	rm $(objs)


