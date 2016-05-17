Copied centralizedCube4GILDAS.lmv-clean from 30Apri16/

To plot RGB CO(2-1) mom0 map:
# make 0th moment seperating out red, blue components Integrated
- total chan range = 126-160, which is consistent with mom0 map now (see 15May16)
blue: channel [160, 144]; velocity [-324.88, 0]
Noise:
    from line-free: 0.53733 Jy km/s/Beam
    calculate: 128.626e-3 Jy km/s/Beam
red: channel [126, 144]; velocity [409, 0]
Noise: 
    from line-free: 0.43225 Jy km/s/Beam
    calculate: 136e-3 Jy km/s/Beam
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RED !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

run map_sum
! in: centralizedCube4GILDAS.lmv-clean
! out: Integrated_red_centralized.lmv-clean
! vel: 0 409      ! z=0 frame

let name Integrated_red_centralized
let type lmv-clean
let size 20
let first 0
let last 0
let spacing 3*0.43225
go nice
pause "Type Continue or Quit"
fits Integrated_red_centralized.fits from Integrated_red_centralized.lmv-clean



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! BLUE !! Map
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
run map_sum
! in: centralizedCube4GILDAS.lmv-clean
! out: Integrated_blue_centralized.lmv-clean
! vel: -651.90 -342.664

let name Integrated_blue_centralized
let type lmv-clean
let size 20
let first 0
let last 0
let spacing 3*128.626e-3
go nice
pause "Type Continue or Quit"
fits Integrated_blue_centralized.fits from Integrated_blue_centralized.lmv-clean