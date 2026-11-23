(* Regenerate the two manuscript figures and their plotted values. *)
Module[{dir = DirectoryName[ExpandFileName[$InputFileName]], raw, headers, rows, take, numeric,
  slice, panels, series, display, makePanel, style, figure1, figure2, write},
  raw = Import[FileNameJoin[{DirectoryName[dir], "Verification", "verification_results.csv"}], "CSV", "Numeric" -> False];
  headers = First[raw]; rows = AssociationThread[headers, #] & /@ Rest[raw];
  take[q_String] := Select[rows, #["q"] == q && #["working_dps"] == "55" &];
  numeric[s_String] := ToExpression[StringReplace[s, "e" -> "*^"]];
  write[name_, cols_, rr_] := Export[FileNameJoin[{dir, name}], Prepend[Lookup[#, cols] & /@ rr, cols], "CSV"];
  slice = take["100"]; panels = take /@ {"100", "400", "1600"};
  If[!AllTrue[panels, Length[#] == 41 && (Rationalize[numeric[#["x"]], 10^-14] & /@ #) === Range[-1, 1, 1/20] &], Return[$Failed]];
  write["figure_1_data.csv", {"q","x","q_abs_R"}, slice];
  write["figure_2_data.csv", {"q","x","abs_one_term","abs_R"}, Join @@ panels];
  style = {FontFamily -> "Times", FontSize -> 12, FontColor -> Black};
  series[rr_, col_] := ({numeric[#["x"]], numeric[#[col]]} & /@ rr);
  figure1 = ListLinePlot[series[slice, "q_abs_R"], Frame -> True, Axes -> False, FrameLabel -> {"x", "q |R|"},
    PlotStyle -> Directive[Black, AbsoluteThickness[1.25]], PlotRange -> {{-1,1},{0,0.086}}, PlotRangePadding -> Scaled[.025],
    PlotLabel -> Style["q = 100", Black], ImageSize -> 410, AspectRatio -> .45, BaseStyle -> style, Background -> White];
  makePanel[rr_] := ListLogPlot[{series[rr, "abs_one_term"],series[rr, "abs_R"]}, Joined -> True,
    PlotStyle -> {Directive[Black,Dashed,AbsoluteThickness[.9]],Directive[Black,AbsoluteThickness[1.1]]},
    PlotMarkers -> {{Style["\[EmptyCircle]",Black],7},{Style["\[FilledSquare]",Black],4}},
    Frame -> True, Axes -> False, PlotRange -> {{-1,1},{0.000005,0.1}},
    FrameTicks -> {{{{0.00001,"0.00001"},{0.0001,"0.0001"},{0.001,"0.001"},{0.01,"0.01"},{0.1,"0.1"}},None},{{-1,0,1},None}},
    FrameLabel -> {"x",None}, PlotLabel -> Style["q = "<>rr[[1]]["q"],Black], ImagePadding -> {{55,9},{38,24}},
    ImageSize -> 236, AspectRatio -> .81, BaseStyle -> style, Background -> White];
  figure2 = Column[{GraphicsRow[makePanel /@ panels, Spacings -> 0, ImageSize -> 724],
    Style["Absolute error    \[EmptyCircle]  one term (dashed)       \[FilledSquare]  two terms (solid)", 12, Black, FontFamily -> "Times"]},
    Alignment -> Center, Background -> White];
  Export[FileNameJoin[{dir,"figure_1.pdf"}], figure1, "PDF"];
  Export[FileNameJoin[{dir,"figure_2.pdf"}], figure2, "PDF"];
  {Length[slice], Total[Length /@ panels]}
]