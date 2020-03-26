# https://medium.com/@porteneuve/mastering-git-subtrees-943d29a798ec

git remote add nim-kmers https://github.com/bio-nim/nim-kmers.git
git remote add nim-networkx https://github.com/bio-nim/nim-networkx.git
git remote add msgpack4nim git@github.com:jangko/msgpack4nim.git
git remote add hts-nim git@github.com:brentp/hts-nim.git
git remote add bio-hts-nim git@github.com:bio-nim/hts-nim.git
git remote add nim-heap git@github.com:bluenote10/nim-heap.git
git remote add nim-heap git@github.com:bio-nim/nim-heap.git
git remote add cligen git@github.com:c-blake/cligen.git
git remote add threadpools git@github.com:yglukhov/threadpools.git
git remote add comprehension git@github.com:alehander92/comprehension.git
git remote add c_alikes git@github.com:ReneSac/c_alikes.git
git remote add BitVector git@github.com:MarcAzar/BitVector.git

git fetch nim-kmers
git fetch nim-networkx

# To start from upstream:
#git read-tree --prefix=vendor/nim-kmers -u nim-kmers/master
#git read-tree --prefix=vendor/nim-networkx -u nim-networkx/nodi

# To update from upstream:
#git merge -s subtree --squash nim-networkx/nodi
## git merge -X subtree=vendor/nim-networkx --squash nim-networkx/nodi
#git commit

# To modify from here:
#git checkout -B backport nim-networkx/nodi
#git cherry-pick -x --strategy=subtree (commits)
# --strategy=subtree is not usually needed.

# To create a new subtree from an existing directory:
#git co -B split
#git filter-branch --subdirectory-filter vendor/mything
#git remote add mything URL
#git push -u mything split:master
