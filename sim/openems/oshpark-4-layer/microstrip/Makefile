.PHONY: post-proc.py
post-proc: post-proc.py sim
	./post-proc.py

sim: proc.py params.py
	./proc.py

.PHONY: csx
csx:
	AppCSXCAD csx.xml

.PHONY: efield
efield:
	paraview
