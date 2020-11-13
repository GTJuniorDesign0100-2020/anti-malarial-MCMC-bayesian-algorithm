import React from 'react';
import JSZip from 'jszip';
import {saveAs} from 'file-saver';

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
          <td>{this.getResultDownloadLinks(csvFileText)}</td>
        </tr>
      )
    })
  }

  getResultDownloadLinks(csvFileText) {
    if (!csvFileText) {
      return '';
    }

    const probability_filename = 'probability_of_recrudescence.csv';

    let extra_filenames = Object.keys(csvFileText);
    // Remove probability file from those included in .zip
    extra_filenames.splice(extra_filenames.indexOf(probability_filename), 1);

    return (
      <div>
        <CSVDownloadLink
          key={probability_filename}
          csvFileName={probability_filename}
          csvFileText={csvFileText[probability_filename]}
        />
        <CSVZipDownloadLink
          fileNames={extra_filenames}
          fileContentDict={csvFileText}
        />
      </div>
    );
  }
}

const CSVZipDownloadLink = ({fileNames, fileContentDict}) => {
  let zip = new JSZip();
  for (let filename of fileNames) {
    const csvText = fileContentDict[filename];
    zip.file(filename, csvText);
  }

  return (
    <p><a
      href="#"
      onClick={(evt) => {
        evt.preventDefault();
        zip.generateAsync({type: 'blob'}).then(
          zipContent => saveAs(zipContent, 'advanced_stats.zip')
        );
      }}>
        advanced_stats.zip
    </a></p>
  );
};

const CSVDownloadLink = ({csvFileName, csvFileText}) => {
  const csvText = `data:text/csv;charset=utf-8,${csvFileText}`;
  const csvUri = encodeURI(csvText);
  return (
    <p><a href={csvUri} download={csvFileName}>{csvFileName}</a></p>
  );
};
