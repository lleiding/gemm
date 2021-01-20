# Hacking GeMM

**Note: this is specifically intended for development of the *Zosterops* project**

*Daniel Vedder, 19/01/2021*

## Important branches

- `master` This should only be updated when there is a new stable (i.e. correctly working) version.

- `zosterops` This is the primary development branch. This doesn't have to be bug-free, but should
  be able to run the model code without crashing.
  
  
## Workflow

1. If you're working on a larger feature, create a feature branch on Github. Otherwise, just use
   the local branch on your own computer.

2. Pull the current development version from [`zosterops`](https://github.com/lleiding/gemm/tree/zosterops).

3. Implement your changes.

4. Run `zrun.sh` to make sure the model executes without crashing.

5. Commit your work frequently.

6. If you're using a local branch, push back to `zosterops` at intervals. Make sure the model runs
   before you do so. Also, do another pull before pushing in case somebody else changed the branch
   in the meantime.
   
7. If you're using a feature branch, push to it as often as you like. When you're done with the
   feature, merge any new developments from `zosterops` into it, then submit a merge request
   back into `zosterops`.
   
8. Repeat :-)

9. Don't break anything important!
