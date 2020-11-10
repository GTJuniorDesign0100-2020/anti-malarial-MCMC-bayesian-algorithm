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
                    <th width='25%'>Date</th>
                    <th width='25%'>Time</th>
                    <th width='25%'>Input File Name</th>
                    <th width='25%'>Output</th>
                  </tr>
                    <tbody>
                        {this.renderTableData(this.props.data)}
                    </tbody>
                </table>
            </div>
        )
    }

    renderTableData(data) {
        return data.map((dataset, index) => {
            const [date, time, fileName, output] = dataset
            return (
                <tr key={date}>
                    <td>{date}</td>
                    <td>{time}</td>
                    <td>{fileName}</td>
                    <td>{output}</td>
                </tr>
            )
        })
    }
}
