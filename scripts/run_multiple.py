import os
import time
import subprocess as sp
import shlex


n = 2500
c = 4
s = 5
nz = 9223372036854775807
ne = 10000

reps = 3

cmds = []
for l in [0.1, 0.2, 0.3, 0.4]:
    for (algo, n_init) in [
        ("bernoulli", 1),
        ("glassermanli", 5),
    ]:
        for rep in range(reps):
            filename = "{}_ninit_{}_n_{}_c_{}_s_{}_l_{}_ne_{}_rep_{}".format(
                algo, n_init, n, c, s, l, ne, rep
            )

            cmd = """
    julia7 src/main.jl --n {n} --c {c} --s {s} --l {l} --nz {nz} --ne {ne} --filename {filename} --n_init {n_init} --a {algo}
            """.format(
                n=n,
                c=c,
                s=s,
                l=l,
                nz=nz,
                ne=ne,
                filename=filename,
                n_init=n_init,
                algo=algo
            )

            cmds.append(cmd)


show = True
if show:
    for cmd in cmds:
        print(cmd)
else:
    procs = []
    for cmd in cmds:
        print(cmd)
        proc = sp.Popen(cmd, shell=True)
        procs.append(proc)

    print("starting to wait for subprocesses...")

    for proc in procs:
        proc.wait()

    print("exiting from main...")
