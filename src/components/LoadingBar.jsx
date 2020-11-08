import React from 'react';
import {Line} from 'rc-progress'

export default function LoadingBar({percentCompleted, secondsRunning}) {
  // Convert from "0 to 1.0" to "0 to 100" scale
  percentCompleted = percentCompleted * 100.0;
  percentCompleted = Math.max(0, Math.min(100, percentCompleted));

  secondsRunning = Math.round(secondsRunning);

  return (
    <div style={{width: "100%", textAlign: "center"}}>
      <p style={{margin: "2px"}}>
        {percentCompleted.toFixed(2)}% (Running for {secondsRunning}s)
      </p>
      <Line
        percent={percentCompleted}
        strokeColor="#84bc49"
        strokeWidth="1"
        trailWidth="1" />
    </div>
  );
}
