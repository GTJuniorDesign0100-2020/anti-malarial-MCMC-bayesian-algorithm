import React, {useState} from 'react';
import LoadingBar from './LoadingBar';

/**
 * TODO: Remove this!
 * A demonstration of how the loading bar component works
 */

export default function DemoLoadingBar() {
    const [startTime, setStartTime] = useState(Date.now());
    const [secondsSinceStart, setSecondsSince] = useState(0.0);

    function updateTime() {
        setSecondsSince((Date.now() - startTime) / 1000.0);
    }
    const percentCompleted = secondsSinceStart / 10.0;

    if (percentCompleted < 1.0) {
        // Update every 10ms
        setTimeout(updateTime, 10);
    }

    const styles = {
        'width': '30vw'
    };

    const content = (percentCompleted < 1.0)
        ? <LoadingBar percentCompleted={percentCompleted} secondsRunning={secondsSinceStart} />
        : 'Completed'

    return (
        <div style={styles}>
            {content}
        </div>
    );
}
