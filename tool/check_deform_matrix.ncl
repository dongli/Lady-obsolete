load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

PI = atan(1.0)*4.0

procedure plot_deform_matrix(wks, H, xs)
begin

    ; plot background
    res_bg = True
    ; res_bg@gsnMaximize = True
    res_bg@gsnDraw = False
    res_bg@gsnFrame = False
    a = 5.0e+06
    res_bg@trXMinF = -a
    res_bg@trXMaxF = a
    res_bg@trYMinF = -a
    res_bg@trYMaxF = a
    res_bg@tmXBTickStartF = -a
    res_bg@tmXBTickEndF = a
    res_bg@tmYLTickStartF = -a
    res_bg@tmYLTickEndF = a

    plot_bg = gsn_blank_plot(wks, res_bg)

    ; plot the skeleton points
    res_xs = True
    res_xs@gsMarkerIndex = 16
    res_xs@gsMarkerColor = "green"
    res_xs@gsMarkerSizeF = 0.005

    plot_xs = gsn_add_polymarker(wks, plot_bg, xs(:,0), xs(:,1), res_xs)

    ; plot deformation matrix
    n = 50
    y = new((/2/), double)
    x = new((/n,2/), double)
    dtheta = 2*PI/(n-1)

    res_x = True
    res_x@gsLineThicknessF = 3

    res_x@gsLineColor = "black"
    do i = 0, n-1
        x(i,:) = H(0,:,:)#(/cos(i*dtheta),sin(i*dtheta)/)
    end do
    plot_x0 = gsn_add_polyline(wks, plot_bg, x(:,0), x(:,1), res_x)

    res_x@gsLineColor = "gray"
    do i = 0, n-1
        x(i,:) = H(1,:,:)#(/cos(i*dtheta),sin(i*dtheta)/)
    end do
    plot_x1 = gsn_add_polyline(wks, plot_bg, x(:,0), x(:,1), res_x)

    res_x@gsLineColor = "gray"
    do i = 0, n-1
        x(i,:) = H(2,:,:)#(/cos(i*dtheta),sin(i*dtheta)/)
    end do
    plot_x2 = gsn_add_polyline(wks, plot_bg, x(:,0), x(:,1), res_x)

    res_x@gsLineColor = "gray"
    do i = 0, n-1
        x(i,:) = H(3,:,:)#(/cos(i*dtheta),sin(i*dtheta)/)
    end do
    plot_x3 = gsn_add_polyline(wks, plot_bg, x(:,0), x(:,1), res_x)

    res_x@gsLineColor = "gray"
    do i = 0, n-1
        x(i,:) = H(4,:,:)#(/cos(i*dtheta),sin(i*dtheta)/)
    end do
    plot_x4 = gsn_add_polyline(wks, plot_bg, x(:,0), x(:,1), res_x)

    draw(plot_bg)
    frame(wks)

end

begin

    load "../build/Debug/deform_matrix.dat"

    wks = gsn_open_wks("pdf", "deform_matrix")

    plot_deform_matrix(wks, H, xs)

end