import React from "react";

export default class Footer extends React.Component {
  constructor(props) {
      super(props);
  }

  render() {
      return (
          <div className="footer">
            <p>To read about the algorithm powering this application, check out its <a target="_blank" rel="noopener noreferrer" href="https://github.com/GTJuniorDesign0100-2020/anti-malarial-MCMC-bayesian-algorithm">article</a> in the American Society for Microbiology.  If you want to contribute to this project or examine the technology, the software is available on our <a target="_blank" rel="noopener noreferrer" href="https://github.com/GTJuniorDesign0100-2020/anti-malarial-MCMC-bayesian-algorithm">GitHub</a>.</p>
          </div>
      );
  }
}
