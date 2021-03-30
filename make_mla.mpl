with(LibraryTools):
loc:= cat(currentdir(),"/MultivariatePowerSeries.mla");
Create(loc);
read "./src/MultivariatePowerSeries.mpl":
eval(MultivariatePowerSeries);
Save(MultivariatePowerSeries, loc);
ShowContents(loc);
