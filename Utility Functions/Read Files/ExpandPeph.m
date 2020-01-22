function Outeph = ExpandPeph(Peph)

Outeph = Peph;

if Peph.NumSV < 32
    prns = Peph.PRN(1:Peph.NumSV);
    
    Outeph.PRN = repmat((1:32)', 1, Peph.NumEpochs);
    Outeph.PRN = reshape(Outeph.PRN, 32*Peph.NumEpochs,1);
    
    Outeph.clock_bias = NaN(32, Peph.NumEpochs);
    tmp = reshape(Peph.clock_bias, Peph.NumSV, Peph.NumEpochs);
    Outeph.clock_bias(prns, :) = tmp;
    Outeph.clock_bias = reshape(Outeph.clock_bias, 32*Peph.NumEpochs,1);

    Outeph.clock_drift = NaN(32, Peph.NumEpochs);
    tmp = reshape(Peph.clock_drift, Peph.NumSV, Peph.NumEpochs);
    Outeph.clock_drift(prns, :) = tmp;
    Outeph.clock_drift = reshape(Outeph.clock_drift, 32*Peph.NumEpochs,1);    
    
    tmppos = NaN(32, Peph.NumEpochs);
    tmp = reshape(Peph.position(:,1), Peph.NumSV, Peph.NumEpochs);
    tmppos(prns, :) = tmp;
    Outeph.position = reshape(tmppos, 32*Peph.NumEpochs,1);    
    
    tmppos = NaN(32, Peph.NumEpochs);
    tmp = reshape(Peph.position(:,2), Peph.NumSV, Peph.NumEpochs);
    tmppos(prns, :) = tmp;
    Outeph.position(:,2) = reshape(tmppos, 32*Peph.NumEpochs,1);    

    tmppos = NaN(32, Peph.NumEpochs);
    tmp = reshape(Peph.position(:,3), Peph.NumSV, Peph.NumEpochs);
    tmppos(prns, :) = tmp;
    Outeph.position(:,3) = reshape(tmppos, 32*Peph.NumEpochs,1);     
    
    tmppos = NaN(32, Peph.NumEpochs);
    tmp = reshape(Peph.velocity(:,1), Peph.NumSV, Peph.NumEpochs);
    tmppos(prns, :) = tmp;
    Outeph.velocity = reshape(tmppos, 32*Peph.NumEpochs,1);    
    
    tmppos = NaN(32, Peph.NumEpochs);
    tmp = reshape(Peph.velocity(:,2), Peph.NumSV, Peph.NumEpochs);
    tmppos(prns, :) = tmp;
    Outeph.velocity(:,2) = reshape(tmppos, 32*Peph.NumEpochs,1);    

    tmppos = NaN(32, Peph.NumEpochs);
    tmp = reshape(Peph.velocity(:,3), Peph.NumSV, Peph.NumEpochs);
    tmppos(prns, :) = tmp;
    Outeph.velocity(:,3) = reshape(tmppos, 32*Peph.NumEpochs,1);     
    
    Outeph.Event = ones(32, Peph.NumEpochs);
    tmp = reshape(Peph.Event, Peph.NumSV, Peph.NumEpochs);
    Outeph.Event(prns, :) = tmp;
    Outeph.Event = reshape(Outeph.Event, 32*Peph.NumEpochs,1);     
    
    Outeph.NumSV = 32;
end
