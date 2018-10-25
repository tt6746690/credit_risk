import os
import time
import subprocess as sp
import shlex
import os

julia = "{}/julia-0.7.0/bin/julia".format(os.environ['HOME'])
# julia = "/Applications/Julia-0.7.app/Contents/Resources/julia/bin/julia"


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
            {julia} src/main.jl --n {n} --c {c} --s {s} --l {l} --nz {nz} --ne {ne} --filename {filename} --n_init {n_init} --a {algo}
            """.format(
                julia=julia,
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


show = False
if show:
    for cmd in cmds:
        print(cmd)
else:
    procs = []
    for cmd in cmds:
        print(cmd)
        proc = sp.Popen(cmd, shell=True, env=os.environ)
        procs.append(proc)

    print("starting to wait for subprocesses...")

    for proc in procs:
        proc.wait()

    print("exiting from main...")
