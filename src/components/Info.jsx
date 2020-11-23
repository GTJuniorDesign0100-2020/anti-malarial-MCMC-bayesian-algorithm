import React from "react";

export default class Info extends React.Component {
  constructor(props) {
      super(props);
  }

  render() {
      return (
          <div className="info">
            <p>To read about the algorithm powering this application, check out its <a target='_blank' href='https://github.com/GTJuniorDesign0100-2020/anti-malarial-MCMC-bayesian-algorithm'>article</a> in the American Society for Microbiology.</p>
          </div>
      );
  }
}
