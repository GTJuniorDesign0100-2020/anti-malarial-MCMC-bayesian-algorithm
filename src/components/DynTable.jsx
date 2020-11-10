import React from 'react';
import LoadingBar from './LoadingBar';
import DemoLoadingBar from './DemoLoadingBar';

export default class DynTable extends React.Component {
    constructor(props) {
        super(props)
    }

    render() {
        return (
            <div>
                <h1 align='center'>Results</h1>
                <table id='results' align='center' border='2px' width='80%'>
                  <tr>
                    <th width='20%'>Date</th>
                    <th width='20%'>Time</th>
                    <th width='20%'>Input File Name</th>
                    <th width='20%'>Status</th>
                    <th width='20%'>Output</th>
                  </tr>
                    <tbody>
                        {this.renderTableData(this.props.data)}
                    </tbody>
                </table>
            </div>
        )
    }

    renderTableData(data) {
        console.log(data);
        return Object.values(data).map((dataset, index) => {
            const {date, inputFilename, loadingBar, results} = dataset;
            console.log(dataset);
            return (
                <tr key={date.toISOString()}>
                    <td>{date.toLocaleDateString()}</td>
                    <td>{date.toLocaleTimeString()}</td>
                    <td>{inputFilename}</td>
                    <td>{loadingBar}</td>
                    <td>TODO</td>
                </tr>
            )
        })
    }
}
