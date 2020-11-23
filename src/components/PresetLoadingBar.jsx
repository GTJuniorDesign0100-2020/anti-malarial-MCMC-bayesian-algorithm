import React, {useState} from 'react';
import LoadingBar from './LoadingBar';

/**
 * A loading bar where the time to complete is known in advance
 */

export default function PresetLoadingBar({timeToComplete=10.0}) {
    const [startTime] = useState(Date.now());
    const [secondsSinceStart, setSecondsSince] = useState(0.0);

    // Set null completion times = 100% completion
    timeToComplete = timeToComplete ? timeToComplete : 0.0001;

    function updateTime() {
        setSecondsSince((Date.now() - startTime) / 1000.0);
    }
    const percentCompleted = secondsSinceStart / timeToComplete;

    // Update every 10ms
    setTimeout(updateTime, 10);

    const styles = {
        'width': '30vw'
    };

    return (
        <div style={styles}>
            <LoadingBar percentCompleted={percentCompleted} secondsRunning={secondsSinceStart} />
        </div>
    );
}
