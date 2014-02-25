load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    start_time = 0
    end_time = 300
    time_step = 300
    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/Debug/tracers.%3.3d.nc\n' $i; done")

    num_time = dimsizes(fs)

    wks = gsn_open_wks("pdf", "cells")

    gsn_define_colormap(wks, "WhBlGrYeRe")

    res = True
    res@cnLinesOn = False
    res@cnFillOn = True
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 0.0
    res@cnMaxLevelValF = 750.0
    res@cnLevelSpacingF = 50.0
    res@gsnPolar = "SH"
    res@mpMaxLatF = -20

    do l = 0, num_time-1
        system("echo time step: "+l)
        f = addfile(fs(l), "r")
        num_lat = dimsizes(f->lat)
        q = f->q1/f->q0
        copy_VarCoords(f->q1, q)
        ;res@gsnRightString = "Max:"+sprintf("%8.1f", max(q(1:num_lat-2,:)))
        plot = gsn_csm_contour_map(wks, q(1:num_lat-2,:), res)
    end do

end
