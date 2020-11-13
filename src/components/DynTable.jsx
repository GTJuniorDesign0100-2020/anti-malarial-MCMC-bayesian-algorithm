import React from 'react';

export default class DynTable extends React.Component {
  constructor(props) {
    super(props)
  }

  render() {
    return (
      <div>
        <h1 align='center'>Results</h1>
        <table id='results' align='center' border='2px' width='80%'>
          <thead>
            <tr>
              <th width='20%'>Date</th>
              <th width='20%'>Time</th>
              <th width='20%'>Input File Name</th>
              <th width='20%'>Status</th>
              <th width='20%'>Output</th>
            </tr>
          </thead>
            <tbody>
              {this.renderTableData(this.props.data)}
            </tbody>
        </table>
      </div>
    )
  }

  renderTableData(data) {
    return Object.values(data).map((dataset) => {
      const {date, inputFilename, status, results} = dataset;
      const csvFileText = results.output_file_text;
      return (
        <tr key={date.toISOString()}>
          <td>{date.toLocaleDateString()}</td>
          <td>{date.toLocaleTimeString()}</td>
          <td>{inputFilename}</td>
          <td>{status}</td>
          <td>{this.getCSVFilesLinks(csvFileText)}</td>
        </tr>
      )
    })
  }

  getCSVFilesLinks(csvFileText) {
    if (!csvFileText) {
      return '';
    }

    return (
      // TODO: Bundle links/files in a zip?
      <div>
        {Object.keys(csvFileText).map(filename =>
          <CSVDownloadLink
            key={filename}
            csvFileName={filename}
            csvFileText={csvFileText[filename]}
          />
        )}
      </div>
    );
  }
}

const CSVDownloadLink = ({csvFileName, csvFileText}) => {
  const csvText = `data:text/csv;charset=utf-8,${csvFileText}`;
  const csvUri = encodeURI(csvText);
  return (
    <p><a href={csvUri} download={csvFileName}>{csvFileName}</a></p>
  );
};
