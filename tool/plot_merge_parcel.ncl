load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    load "./merge_parcel.txt"
    load "/Users/dongli/Dropbox/Works/Lady/build/Debug/merge_parcel.txt"

    PI = atan(1.0)*4
    RAD = PI/180
    point_size = 0.1
    line_thickness = 3
    satellite_angle = 3

    wks = gsn_open_wks("pdf", "merge_parcel")

    res_map = True
    res_map@gsnMaximize = False
    res_map@gsnFrame = False
    res_map@mpCenterLonF = 180.0
    ; res_map@mpProjection = "Satellite"
    ; res_map@mpCenterLonF = centroid0(0)
    ; res_map@mpCenterLatF = centroid0(1)

    ; res_map@mpLimitMode = "Angles"
    ; res_map@mpLeftAngleF = satellite_angle
    ; res_map@mpRightAngleF = satellite_angle
    ; res_map@mpTopAngleF = satellite_angle
    ; res_map@mpBottomAngleF = satellite_angle
    res_map@mpLimitMode = "LatLon"
    res_map@mpMinLonF = 240
    res_map@mpMaxLonF = 270
    res_map@mpMinLatF = -55
    res_map@mpMaxLatF = -25
    res_map@pmTickMarkDisplayMode = "Always"

    res_map@mpFillOn = True
    res_map@mpOutlineOn = True
    res_map@mpGridAndLimbOn = True
    res_map@mpGridLineColor = "gray"
    res_map@mpGridLonSpacingF = 1.5
    res_map@mpGridLatSpacingF = 1.5
    map = gsn_csm_map(wks, res_map)

    res_shape = True
    res_shape@gsMarkerIndex = 16
    res_shape@gsMarkerSizeF = point_size
    res_shape@gsLineThicknessF = line_thickness

    gsn_define_colormap(wks, "MPL_Reds")
    ; gsn_reverse_colormap(wks)
    colormap = gsn_retrieve_colormap(wks)
    levels = fspan(0, 1, 101)
    colors = span_color_indexes(colormap(3:120,:), dimsizes(levels))
print(weights)
    if (isdefined("shape1")) then
        res_shape@gsLineColor = colors(closest_val(weights(0), levels))+3
        gsn_polyline(wks, map, shape1(:,0), shape1(:,1), res_shape)
    end if
    if (isdefined("shape2")) then
        res_shape@gsLineColor = colors(closest_val(weights(1), levels))+3
        gsn_polyline(wks, map, shape2(:,0), shape2(:,1), res_shape)
    end if
    if (isdefined("shape3")) then
        res_shape@gsLineColor = colors(closest_val(weights(2), levels))+3
        gsn_polyline(wks, map, shape3(:,0), shape3(:,1), res_shape)
    end if
    if (isdefined("shape4")) then
        res_shape@gsLineColor = colors(closest_val(weights(3), levels))+3
        gsn_polyline(wks, map, shape4(:,0), shape4(:,1), res_shape)
    end if
    if (isdefined("shape5")) then
        res_shape@gsLineColor = colors(closest_val(weights(4), levels))+3
        gsn_polyline(wks, map, shape5(:,0), shape5(:,1), res_shape)
    end if
    if (isdefined("shape7")) then
        res_shape@gsLineColor = colors(closest_val(weights(6), levels))+3
        gsn_polyline(wks, map, shape7(:,0), shape7(:,1), res_shape)
    end if
    if (isdefined("shape8")) then
        res_shape@gsLineColor = colors(closest_val(weights(7), levels))+3
        gsn_polyline(wks, map, shape8(:,0), shape8(:,1), res_shape)
    end if
    if (isdefined("shape9")) then
        res_shape@gsLineColor = colors(closest_val(weights(8), levels))+3
        gsn_polyline(wks, map, shape9(:,0), shape9(:,1), res_shape)
    end if
    if (isdefined("shape10")) then
        res_shape@gsLineColor = colors(closest_val(weights(9), levels))+3
        gsn_polyline(wks, map, shape10(:,0), shape10(:,1), res_shape)
    end if
    if (isdefined("shape11")) then
        res_shape@gsLineColor = colors(closest_val(weights(10), levels))+3
        gsn_polyline(wks, map, shape11(:,0), shape11(:,1), res_shape)
    end if
    if (isdefined("shape12")) then
        res_shape@gsLineColor = colors(closest_val(weights(11), levels))+3
        gsn_polyline(wks, map, shape12(:,0), shape12(:,1), res_shape)
    end if
    if (isdefined("shape13")) then
        res_shape@gsLineColor = colors(closest_val(weights(12), levels))+3
        gsn_polyline(wks, map, shape13(:,0), shape13(:,1), res_shape)
    end if
    if (isdefined("shape14")) then
        res_shape@gsLineColor = colors(closest_val(weights(13), levels))+3
        gsn_polyline(wks, map, shape14(:,0), shape14(:,1), res_shape)
    end if
    if (isdefined("shape15")) then
        res_shape@gsLineColor = colors(closest_val(weights(14), levels))+3
        gsn_polyline(wks, map, shape15(:,0), shape15(:,1), res_shape)
    end if
    if (isdefined("shape16")) then
        res_shape@gsLineColor = colors(closest_val(weights(15), levels))+3
        gsn_polyline(wks, map, shape16(:,0), shape16(:,1), res_shape)
    end if
    if (isdefined("shape17")) then
        res_shape@gsLineColor = colors(closest_val(weights(16), levels))+3
        gsn_polyline(wks, map, shape17(:,0), shape17(:,1), res_shape)
    end if
    if (isdefined("shape18")) then
        res_shape@gsLineColor = colors(closest_val(weights(17), levels))+3
        gsn_polyline(wks, map, shape18(:,0), shape18(:,1), res_shape)
    end if
    if (isdefined("shape19")) then
        res_shape@gsLineColor = colors(closest_val(weights(18), levels))+3
        gsn_polyline(wks, map, shape19(:,0), shape19(:,1), res_shape)
    end if


    res_shape@gsLineColor = colors(closest_val(weights2(0), levels))+3
    gsn_polyline(wks, map, shape101(:,0), shape101(:,1), res_shape)
    res_shape@gsLineColor = colors(closest_val(weights2(1), levels))+3
    gsn_polyline(wks, map, shape102(:,0), shape102(:,1), res_shape)
    res_shape@gsLineColor = colors(closest_val(weights2(2), levels))+3
    gsn_polyline(wks, map, shape103(:,0), shape103(:,1), res_shape)
    res_shape@gsLineColor = colors(closest_val(weights2(3), levels))+3
    gsn_polyline(wks, map, shape104(:,0), shape104(:,1), res_shape)
    res_shape@gsLineColor = colors(closest_val(weights2(4), levels))+3
    gsn_polyline(wks, map, shape105(:,0), shape105(:,1), res_shape)
    res_shape@gsLineColor = colors(closest_val(weights2(5), levels))+3
    gsn_polyline(wks, map, shape106(:,0), shape106(:,1), res_shape)
    res_shape@gsLineColor = colors(closest_val(weights2(6), levels))+3
    gsn_polyline(wks, map, shape107(:,0), shape107(:,1), res_shape)
    res_shape@gsLineColor = colors(closest_val(weights2(7), levels))+3
    gsn_polyline(wks, map, shape108(:,0), shape108(:,1), res_shape)
    res_shape@gsLineColor = colors(closest_val(weights2(9), levels))+3
    gsn_polyline(wks, map, shape110(:,0), shape110(:,1), res_shape)
    if (isdefined("shape6")) then
        res_shape@gsLineColor = colors(closest_val(weights(5), levels))+3
        gsn_polyline(wks, map, shape6(:,0), shape6(:,1), res_shape)
    end if
    res_shape@gsLineColor = colors(closest_val(weights2(8), levels))+3
    gsn_polyline(wks, map, shape109(:,0), shape109(:,1), res_shape)


    delete(res_shape@gsLineColor)
    res_shape@gsLineColor = "gray"
    gsn_polyline(wks, map, shape100(:,0), shape100(:,1), res_shape)
    gsn_polyline(wks, map, shape0(:,0), shape0(:,1), res_shape)

    frame(wks)

end