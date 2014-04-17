load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    start_time = 0
    end_time = 101
    time_step = 1
    file_prefix = "tracers.barotropic.80x40"

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/Debug/"+file_prefix+".%3.3d.nc\n' $i; done")

    num_time = dimsizes(fs)

    wks = gsn_open_wks("pdf", file_prefix+".cells")

    ; gsn_define_colormap(wks, "WhBlGrYeRe")
    ; gsn_define_colormap(wks, "ViBlGrWhYeOrRe")
    ; gsn_define_colormap(wks, "BlueYellowRed")
    cmap = (/(/1.0, 1.0, 1.0 /),(/0.0, 0.0, 0.0 /), \
             (/0.84,0.84,1.0 /),(/0.38,0.4, 0.98/), \
             (/0.04,0.15,0.89/),(/0.06,0.4, 0.46/), \
             (/0.15,0.71,0.14/),(/0.54,0.85,0.17/), \
             (/0.98,0.98,0.22/),(/0.99,0.53,0.14/), \
             (/0.99,0.09,0.11/)/)
    gsn_define_colormap(wks, cmap)

    res = True
    ;res@gsnMaximize = True
    res@cnLinesOn = False
    res@cnFillOn = True
    ; res@cnFillMode = "RasterFill"
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ManualLevels"
    ; res@cnMinLevelValF = 75000
    ; res@cnMaxLevelValF = 85000
    ; res@cnLevelSpacingF = 500.0
    res@cnMinLevelValF = 0.2
    ; res@cnMaxLevelValF = 1000.0
    ; res@cnLevelSpacingF = 50.0
    res@cnMaxLevelValF = 0.9
    res@cnLevelSpacingF = 0.1
    ; res@cnMinLevelValF = 0
    ; res@cnMaxLevelValF = 2
    ; res@cnLevelSpacingF = 0.1
    ; res@gsnPolar = "NH"
    ; res@mpMinLatF = 60
    res@mpCenterLonF = 180.0
    res@mpOutlineOn = False
    res@pmTickMarkDisplayMode = "Always"

    res@lbBoxMinorExtentF = 0.2

    do l = 0, num_time-1
        system("echo file: "+fs(l))
        f = addfile(fs(l), "r")
        num_lat = dimsizes(f->lat)
        ; q = f->q0
        q = f->q1/f->q0
        copy_VarCoords(f->q0, q)
        ;res@gsnRightString = "Max:"+sprintf("%8.1f", max(q(1:num_lat-2,:)))
        ;plot = gsn_csm_contour_map(wks, q(1:num_lat-2,:), res)
        res@gsnLeftString = ""
        res@gsnRightString = "Min:"+sprintf("%8.5f", min(q))+" "+ \
                             "Max:"+sprintf("%8.5f", max(q))
        plot = gsn_csm_contour_map(wks, q(1:num_lat-2,:), res)
    end do

end
