function [delay] = saastamoinen_model_SU(lat, lon, h, el)
%         global geoid

        %Saastamoinen model requires (positive) orthometric height
%         if (exist('geoid','var') && isfield(geoid,'ncols') && geoid.ncols ~= 0)
%             %geoid ondulation interpolation
%             undu = grid_bilin_interp(lon, lat, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
%             h = h - undu;
%         end
        h(h < 0) = 0;

        if (h < 5000)

            %conversion to radians
            el = abs(el) * pi/180;

            %Standard atmosphere - Berg, 1948 (Bernese)
            %pressure [mbar]
            Pr =  1013.25;% goGNSS.STD_PRES;
            %temperature [K]
            Tr = 291.15;%goGNSS.STD_TEMP;
            %humidity [%]
            Hr = 50;%goGNSS.STD_HUMI;

            P = Pr * (1-0.0000226*h).^5.225;
            T = Tr - 0.0065*h;
            H = Hr * exp(-0.0006396*h);

            %----------------------------------------------------------------------

            %linear interpolation
            h_a = [0; 500; 1000; 1500; 2000; 2500; 3000; 4000; 5000];
            B_a = [1.156; 1.079; 1.006; 0.938; 0.874; 0.813; 0.757; 0.654; 0.563];

            t = zeros(length(T),1);
            B = zeros(length(T),1);

            for i = 1 : length(T)

                d = h_a - h(i);
                [~, j] = min(abs(d));
                if (d(j) > 0)
                    index = [j-1; j];
                else
                    index = [j; j+1];
                end

                t(i) = (h(i) - h_a(index(1))) ./ (h_a(index(2)) - h_a(index(1)));
                B(i) = (1-t(i))*B_a(index(1)) + t(i)*B_a(index(2));
            end

            %----------------------------------------------------------------------

            e = 0.01 * H .* exp(-37.2465 + 0.213166*T - 0.000256908*T.^2);

            %tropospheric delay
            delay = ((0.002277 ./ sin(el)) .* (P - (B ./ (tan(el)).^2)) + (0.002277 ./ sin(el)) .* (1255./T + 0.05) .* e);
        else
            delay = zeros(size(el));
        end
    end