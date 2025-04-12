function anovaTbl = bbar_partialEta(anovaTbl)

    for i = 1:height(anovaTbl)
    F = anovaTbl.FStat(i);
    df1 = anovaTbl.DF1(i);
    df2 = anovaTbl.DF2(i);
    pes = (F * df1) / (F * df1 + df2);
    anovaTbl.PartialEtaSq(i) = pes;
    % fprintf('Effect %d: Partial Eta-squared = %.4f\n', i, pes);
    end

end