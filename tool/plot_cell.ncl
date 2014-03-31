load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    start_time = 0
    end_time = 0
    time_step = 1
    file_prefix = "tracers.barotropic.80x40"

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/Debug/"+file_prefix+".%3.3d.nc\n' $i; done")

    num_time = dimsizes(fs)

    wks = gsn_open_wks("pdf", file_prefix+".cells")

    gsn_define_colormap(wks, "WhBlGrYeRe")
    ; gsn_define_colormap(wks, "ViBlGrWhYeOrRe")

    res = True
    ;res@gsnMaximize = True
    res@cnLinesOn = False
    res@cnFillOn = True
    res@cnFillMode = "RasterFill"
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 78000
    res@cnMaxLevelValF = 88000
    res@cnLevelSpacingF = 500.0
    ; res@cnMaxLevelValF = 1000.0
    ; res@cnLevelSpacingF = 50.0
    ; res@cnMaxLevelValF = 1.1
    ; res@cnLevelSpacingF = 0.05
    ; res@cnMinLevelValF = 0
    ; res@cnMaxLevelValF = 2
    ; res@cnLevelSpacingF = 0.1
    ; res@gsnPolar = "NH"
    ; res@mpMinLatF = 30
    res@mpCenterLonF = 180.0

    do l = 0, num_time-1
        system("echo time step: "+l)
        f = addfile(fs(l), "r")
        num_lat = dimsizes(f->lat)
        q = f->q0
        ; q = f->q1/f->q0
        copy_VarCoords(f->q0, q)
        ;res@gsnRightString = "Max:"+sprintf("%8.1f", max(q(1:num_lat-2,:)))
        ;plot = gsn_csm_contour_map(wks, q(1:num_lat-2,:), res)
        res@gsnRightString = "Max:"+sprintf("%8.1f", max(q))
        plot = gsn_csm_contour_map(wks, q, res)
    end do

end
