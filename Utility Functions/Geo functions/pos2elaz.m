function [el,az] = pos2elaz(usrPos,satPos)

llh = xyz2llh(usrPos);

RxyzEnu = findxyz2enu(llh(1)*pi/180,llh(2)*pi/180);
usr_ehat = RxyzEnu(:,1);
usr_nhat = RxyzEnu(:,2);
usr_uhat = RxyzEnu(:,3);


losxyzb = find_los_xyzb(usrPos,satPos);

los_enub=calc_los_enub(losxyzb,usr_ehat',usr_nhat',usr_uhat');
[el, az] = find_elaz(los_enub);



end