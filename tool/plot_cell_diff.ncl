load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    run_type = "Debug"
    time_stamp = "%4.4d"
    start_time = 150
    end_time = 150
    time_step = 1
    file_prefix_1 = "tracers.deform.sc.240x120.split100"
    file_prefix_2 = "tracers.deform.sc.240x120.split20.2"

    fs_1 = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/"+run_type+"/"+file_prefix_1+"."+time_stamp+".nc\n' $i; done")
    fs_2 = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/"+run_type+"/"+file_prefix_2+"."+time_stamp+".nc\n' $i; done")

    num_time = dimsizes(fs_1)

    wks = gsn_open_wks("pdf", file_prefix_1+"_vs_"+file_prefix_2+".cells")

    ; gsn_define_colormap(wks, "WhBlGrYeRe")
    gsn_define_colormap(wks, "ViBlGrWhYeOrRe")

    res = True
    res@cnLinesOn = False
    res@cnFillOn = True
    res@cnFillMode = "RasterFill"
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -0.01
    res@cnMaxLevelValF = 0.01
    res@cnLevelSpacingF = 0.001
    res@mpCenterLonF = 180.0
    res@mpOutlineOn = False
    res@pmTickMarkDisplayMode = "Always"

    do l = 0, num_time-1
        system("echo file: "+fs_1(l)+" vs "+fs_2(l))
        f_1 = addfile(fs_1(l), "r")
        f_2 = addfile(fs_2(l), "r")
        num_lat = dimsizes(f_1->lat)
        q1 = f_1->q1
        q2 = f_2->q1
        dq = q2-q1
        copy_VarCoords(f_1->q0, dq)
        res@gsnRightString = "Max:"+sprintf("%10.8f", max(dq))
        plot = gsn_csm_contour_map(wks, dq, res)
    end do

end
