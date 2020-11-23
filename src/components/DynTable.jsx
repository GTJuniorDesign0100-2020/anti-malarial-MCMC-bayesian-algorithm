import React from 'react';
import SortableTable from './SortableTable';

export default function DynTable(props) {
  function renderResultRow(resultData) {
    const {date, inputFilename, status, results} = resultData;
    const csvFileText = results.output_file_text;
    return (
      <tr key={date.toISOString()}>
        <td>{date.toLocaleDateString()}</td>
        <td>{date.toLocaleTimeString()}</td>
        <td>{inputFilename}</td>
        <td>{status}</td>
        <td>{getCSVFilesLinks(csvFileText)}</td>
      </tr>
    );
  }

  return (
    <div>
      <h2 align='center'>Results</h2>

      <SortableTable
        columnNames={['Date', 'Time', 'Input File Name', 'Status', 'Output']}
        columnSortKeys={['date', 'date', 'inputFilename', 'status', '']}
        items={
          // Reverse so newest results appear at the top
          Object.values(props.data).reverse()
        }
        itemToTableRowFunc={renderResultRow}
      />
    </div>
  );
}

function getCSVFilesLinks(csvFileText) {
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

const CSVDownloadLink = ({csvFileName, csvFileText}) => {
  const csvText = `data:text/csv;charset=utf-8,${csvFileText}`;
  const csvUri = encodeURI(csvText);
  return (
    <p><a href={csvUri} download={csvFileName}>{csvFileName}</a></p>
  );
};
