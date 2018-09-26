RSYNC = /usr/local/Cellar/rsync/3.1.3_1/bin/rsync
RSYNCTAGS = --archive --verbose --info=progress2 -au
SRC_FOLDER = $(HOME)/github/credit_risk

REMOTE = wpq@comps0.cs.toronto.edu

synccode:
	$(RSYNC) $(RSYNCTAGS) $(SRC_FOLDER) $(REMOTE):/u/wpq/github/credit_risk

data = /h/96/wpq/github/credit_risk/*.txt
syncdata:
	$(RSYNC) $(RSYNCTAGS) $(REMOTE):$(data) $(SRC_FOLDER)/

synccoderev:
	$(RSYNC) $(RSYNCTAGS) $(REMOTE):/u/wpq/github/credit_risk/ $(HOME)/github/credit_risk

avoidpasswordduringssh:
	$(RSYNC) $(RSYNCTAGS) ~/.ssh/id_rsa.pub $(REMOTE):/u/wpq/.ssh

clean:
	rm src/*.mem ||:
	# rm *.pdf     ||:
