# Notebook

## 11 July 2022

Videos for 

## 30 May 2022

Candidate video for frame grab for paper:

- open field, 20190613_1054_icab_hdr_R

Error found here due to HNI moving to the right before the labels start. Now corrected, but will need to run entire process again:

- novel_object/20190611_1552_icab_hni_R

Also noticed this video is too short.

- open_field/20190612_1236_icab_kaga_R_of

Removed all split videos and need to re-run.

## 23 May 2022

Noticed all the fishes were incorrectly labelled, and corrected.

Double-check:

- 20190611_1552_icab_hni_R NO
- 20190616_1227_icab_kaga_R NO (q1)
- 20190616_1506_icab_hni_L NO (q1)

## 4 May 2022

Discussion with Ewan re: paper:

* Time + spatial dependence: plot time dependance of iCab control vs iCab/HdrR -- states 1:3?

* How to show we tracked fish properly? Check random sample. 

* Time dependence -- calculate densities so that they sum to 1 for every second

* Choose subset of states re: time and spatial dependence, then DGE

* Variance explained: inverse-normalise proportion of time spent in each state, then use additive model (e.g. state ~ quadrant + time). Run separately on assays.

* Discussion: 

    - Read Hideaki's papers

    - Cite Cornelius's HMM paper

    - Say going to run over the MIKK panel

## 29 April 2022

Ewan's ideas for figures:

* Figure 1 diagram of the set up, the pairing and then some tracks

* Figure 2 is something about “our tracking is good”

* Figure 3 polar plots and boxplots for DGE, plus HMM tiles

* Figure 4 is the same for SGE

* One thought on the tracks, there is a big time dependence of some of the states (most obviously the still states at the start) - I was wondering if we should pick a couple of states which show strongest time dependence (and / or time by strain dependence)

* My other thought is whether there is a spatial as well as time dependence of the HMM states - eg, are some states the “close the wall” state and others in the middle

* These two things go, along with the genetics, to “the states are meaningful”

## 16 March 2022

Tracking success:

* 188/614 (30%) perfectly tracked

* 483/614 (78%) > 99% tracked

* Minimum is 85%

Progress to data analysis.

## 20220311

Tracking success: 

* 120/614 (20%) perfectly tracked

* 383/614 (62%) > 99% tracked

Most tracking errors caused by lack of background subtraction.
Make bgsub = True as the default, as these inbred lines tend to move around a lot relative to MIKK, which makes the "medaka ghost" phenomenon when bgsub = True less likely.

