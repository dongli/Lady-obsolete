load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    run_type = "Release"
    test_case = "deform"
    ic = ".sc"
    time_stamp = "%4.4d"
    start_time = 300
    end_time = 300
    time_step = 1
    file_prefix_1 = "tracers."+test_case+".240x120.reset_skel.split100"
    file_prefix_2 = "tracers."+test_case+".240x120.no_merge.split20"

    fs1 = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/"+run_type+"/"+file_prefix_1+"."+time_stamp+".nc\n' $i; done")
    fs2 = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/"+run_type+"/"+file_prefix_2+"."+time_stamp+".nc\n' $i; done")

    num_time = dimsizes(fs1)

    wks = gsn_open_wks("pdf", file_prefix_1+"_vs_"+file_prefix_2+ic+".cells")

    gsn_define_colormap(wks, "ViBlGrWhYeOrRe")

    res = True
    res@cnLinesOn = False
    res@cnFillOn = True
    res@cnFillMode = "RasterFill"
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -0.1
    res@cnMaxLevelValF = 0.1
    res@cnLevelSpacingF = 0.01
    res@mpCenterLonF = 180.0
    res@mpOutlineOn = False
    res@pmTickMarkDisplayMode = "Always"
    res@lbBoxMinorExtentF = 0.2

    do l = 0, num_time-1
        system("echo file: "+fs1(l)+" vs "+fs2(l))
        f1 = addfile(fs1(l), "r")
        f2 = addfile(fs2(l), "r")
        if (isStrSubset(test_case, "rotation")) then
            q1 = f1->q1/f1->q0
            q2 = f2->q1/f2->q0
        end if
        if (test_case .eq. "deform") then
            if (ic .eq. ".ch") then
                q1 = f1->q1/f1->q0
                q2 = f2->q1/f2->q0
            end if
            if (ic .eq. ".sc") then
                q1 = f1->q3/f1->q0
                q2 = f2->q3/f2->q0
            end if
            if (ic .eq. ".gh") then
                q1 = f1->q4/f1->q0
                q2 = f2->q4/f2->q0
            end if
        end if
        dq = q2-q1
        copy_VarCoords(f1->q0, dq)
        res@gsnRightString = "Max:"+sprintf("%6.3f", max(dq))
        plot = gsn_csm_contour_map(wks, dq, res)
    end do

end
